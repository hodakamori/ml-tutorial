"""
SOAP Client Demo Application
SOAPサーバーに接続して各種サービスを呼び出すクライアント
"""

from zeep import Client
from zeep.exceptions import Fault


def test_calculator_service(client):
    """計算サービスのテスト"""
    print("=" * 60)
    print("Calculator Service Test")
    print("=" * 60)

    try:
        # 足し算
        result = client.service.add(10, 5)
        print(f"add(10, 5) = {result}")

        # 引き算
        result = client.service.subtract(10, 5)
        print(f"subtract(10, 5) = {result}")

        # 掛け算
        result = client.service.multiply(10, 5)
        print(f"multiply(10, 5) = {result}")

        # 割り算
        result = client.service.divide(10.0, 5.0)
        print(f"divide(10.0, 5.0) = {result}")

        # ゼロ除算のテスト
        try:
            result = client.service.divide(10.0, 0.0)
            print(f"divide(10.0, 0.0) = {result}")
        except Fault as e:
            print(f"divide(10.0, 0.0) -> Error: {e.message}")

    except Exception as e:
        print(f"Error in calculator service: {e}")


def test_user_service(client):
    """ユーザー情報サービスのテスト"""
    print("\n" + "=" * 60)
    print("User Service Test")
    print("=" * 60)

    try:
        # ユーザー1の情報を取得
        user_id = 1
        name = client.service.get_user_name(user_id)
        email = client.service.get_user_email(user_id)
        age = client.service.get_user_age(user_id)
        print(f"User {user_id}: {name}, {email}, {age}歳")

        # ユーザー2の情報を取得
        user_id = 2
        name = client.service.get_user_name(user_id)
        email = client.service.get_user_email(user_id)
        age = client.service.get_user_age(user_id)
        print(f"User {user_id}: {name}, {email}, {age}歳")

        # ユーザー3の情報を取得
        user_id = 3
        name = client.service.get_user_name(user_id)
        email = client.service.get_user_email(user_id)
        age = client.service.get_user_age(user_id)
        print(f"User {user_id}: {name}, {email}, {age}歳")

        # 存在しないユーザー
        user_id = 999
        name = client.service.get_user_name(user_id)
        print(f"User {user_id}: {name}")

    except Exception as e:
        print(f"Error in user service: {e}")


def main():
    """メイン関数"""
    # SOAPサーバーのWSDL URL
    WSDL_URL = 'http://localhost:8000/?wsdl'

    print("Connecting to SOAP server...")
    print(f"WSDL URL: {WSDL_URL}\n")

    try:
        # SOAPクライアントを作成
        client = Client(WSDL_URL)

        # 利用可能なサービスを表示
        print("Available services:")
        for service in client.wsdl.services.values():
            print(f"  - {service.name}")
            for port in service.ports.values():
                for operation in port.binding._operations.values():
                    print(f"    * {operation.name}")
        print()

        # 各サービスをテスト
        test_calculator_service(client)
        test_user_service(client)

        print("\n" + "=" * 60)
        print("All tests completed!")
        print("=" * 60)

    except Exception as e:
        print(f"Error: {e}")
        print("\nMake sure the SOAP server is running!")
        print("Start the server with: python soap_server.py")


if __name__ == '__main__':
    main()
